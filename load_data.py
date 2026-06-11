#!/usr/bin/env python3
"""Initialize a SQLite database from cell-count.csv and load all rows.

Run directly with no arguments:

    python load_data.py

This creates ``cell-count.db`` in the repository root with a normalized
schema and loads every row from ``cell-count.csv``.

Schema (3 normalized tables)
----------------------------
subjects     One row per patient. Holds the attributes that are constant
             across all of a subject's samples (project, condition, age,
             sex, treatment, response).
samples      One row per biological sample, linked to a subject. Holds the
             sample-level attributes (sample_type, time_from_treatment_start).
cell_counts  Long/tidy format: one row per (sample, population) pair holding
             the raw count. Storing cell populations as rows rather than five
             fixed columns lets downstream analysis compute relative
             frequencies generically and makes adding populations trivial.
"""

import csv
import os
import sqlite3

CSV_PATH = "cell-count.csv"
DB_PATH = "cell-count.db"

# The five immune cell populations stored long-format in cell_counts.
POPULATIONS = ["b_cell", "cd8_t_cell", "cd4_t_cell", "nk_cell", "monocyte"]

SCHEMA = """
CREATE TABLE subjects (
    subject   TEXT PRIMARY KEY,
    project   TEXT NOT NULL,
    condition TEXT NOT NULL,        -- indication: melanoma, carcinoma, healthy
    age       INTEGER,
    sex       TEXT,                 -- M / F
    treatment TEXT,                 -- miraclib, phauximab, none
    response  TEXT                  -- yes / no / NULL (e.g. untreated healthy)
);

CREATE TABLE samples (
    sample                    TEXT PRIMARY KEY,
    subject                   TEXT NOT NULL,
    sample_type               TEXT,        -- PBMC / WB
    time_from_treatment_start INTEGER,     -- days: 0, 7, 14
    FOREIGN KEY (subject) REFERENCES subjects (subject)
);

CREATE TABLE cell_counts (
    sample     TEXT NOT NULL,
    population TEXT NOT NULL,       -- b_cell, cd8_t_cell, cd4_t_cell, nk_cell, monocyte
    count      INTEGER NOT NULL,
    PRIMARY KEY (sample, population),
    FOREIGN KEY (sample) REFERENCES samples (sample)
);

-- Indexes to speed up the typical filter/group-by queries used downstream.
CREATE INDEX idx_samples_subject     ON samples (subject);
CREATE INDEX idx_subjects_condition  ON subjects (condition);
CREATE INDEX idx_subjects_treatment  ON subjects (treatment);
CREATE INDEX idx_cell_counts_pop     ON cell_counts (population);
"""


def _blank_to_none(value):
    """Treat empty / whitespace-only CSV cells as SQL NULL."""
    if value is None:
        return None
    value = value.strip()
    return value if value else None


def _to_int(value):
    value = _blank_to_none(value)
    return int(value) if value is not None else None


def load(csv_path=CSV_PATH, db_path=DB_PATH):
    if not os.path.exists(csv_path):
        raise FileNotFoundError(f"Cannot find {csv_path} in the current directory.")

    # Start from a clean database so the script is safely re-runnable.
    if os.path.exists(db_path):
        os.remove(db_path)

    conn = sqlite3.connect(db_path)
    try:
        conn.execute("PRAGMA foreign_keys = ON")
        conn.executescript(SCHEMA)

        subjects = {}   # subject -> row tuple (deduped)
        samples = []    # sample row tuples
        cell_rows = []  # (sample, population, count) tuples

        with open(csv_path, newline="") as f:
            reader = csv.DictReader(f)
            for row in reader:
                subj = row["subject"]
                if subj not in subjects:
                    subjects[subj] = (
                        subj,
                        row["project"],
                        row["condition"],
                        _to_int(row["age"]),
                        _blank_to_none(row["sex"]),
                        _blank_to_none(row["treatment"]),
                        _blank_to_none(row["response"]),
                    )

                sample_id = row["sample"]
                samples.append((
                    sample_id,
                    subj,
                    _blank_to_none(row["sample_type"]),
                    _to_int(row["time_from_treatment_start"]),
                ))

                for pop in POPULATIONS:
                    cell_rows.append((sample_id, pop, _to_int(row[pop])))

        conn.executemany(
            "INSERT INTO subjects VALUES (?, ?, ?, ?, ?, ?, ?)",
            subjects.values(),
        )
        conn.executemany(
            "INSERT INTO samples VALUES (?, ?, ?, ?)",
            samples,
        )
        conn.executemany(
            "INSERT INTO cell_counts VALUES (?, ?, ?)",
            cell_rows,
        )
        conn.commit()

        print(f"Database created at: {db_path}")
        print(f"  subjects:    {len(subjects):,}")
        print(f"  samples:     {len(samples):,}")
        print(f"  cell_counts: {len(cell_rows):,}")
    finally:
        conn.close()


if __name__ == "__main__":
    load()
