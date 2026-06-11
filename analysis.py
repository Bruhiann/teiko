#!/usr/bin/env python3
"""Analytics layer for the cell-count database.

This module holds the reusable analysis functions that answer Bob's
questions. Each function takes an open SQLite connection (or opens the
default ``cell-count.db``) and returns a pandas DataFrame, so the same
logic backs both the command line and the interactive dashboard.

Run directly to print/save the Part 2 cell-frequency summary table:

    python analysis.py
"""

import sqlite3

import pandas as pd

DB_PATH = "cell-count.db"


def connect(db_path=DB_PATH):
    """Open a read-only-style connection to the database."""
    return sqlite3.connect(db_path)


# ---------------------------------------------------------------------------
# Part 2 - Data overview: relative frequency of each population per sample
# ---------------------------------------------------------------------------
def cell_frequencies(conn=None, db_path=DB_PATH):
    """Return per-sample relative frequency of each cell population.

    One row per (sample, population) with columns:
        sample       - sample id
        total_count  - total cells in the sample (sum across all populations)
        population   - immune cell population name
        count        - raw count for that population
        percentage   - count as a percentage of total_count

    The per-sample total is computed with a window function so every
    population row carries its sample's total, then percentage is
    count / total * 100.
    """
    close = conn is None
    conn = conn or connect(db_path)
    try:
        query = """
            SELECT
                sample,
                SUM(count) OVER (PARTITION BY sample)                 AS total_count,
                population,
                count,
                100.0 * count / SUM(count) OVER (PARTITION BY sample) AS percentage
            FROM cell_counts
            ORDER BY sample, population
        """
        df = pd.read_sql_query(query, conn)
    finally:
        if close:
            conn.close()
    return df


if __name__ == "__main__":
    freq = cell_frequencies()

    # Console preview with percentage rounded for readability.
    preview = freq.copy()
    preview["percentage"] = preview["percentage"].round(2)
    print(f"Cell frequency summary table: {len(freq):,} rows "
          f"({freq['sample'].nunique():,} samples x "
          f"{freq['population'].nunique()} populations)\n")
    print(preview.head(15).to_string(index=False))
    print("...")

    # Sanity check: each sample's percentages must sum to 100.
    sums = freq.groupby("sample")["percentage"].sum()
    print(f"\nPer-sample percentage sums: min={sums.min():.4f}, "
          f"max={sums.max():.4f} (expected 100.0)")

    out = "cell_frequencies.csv"
    freq.to_csv(out, index=False)
    print(f"Full table written to {out}")
