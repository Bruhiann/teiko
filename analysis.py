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
from scipy.stats import mannwhitneyu, false_discovery_control

DB_PATH = "cell-count.db"

POPULATIONS = ["b_cell", "cd8_t_cell", "cd4_t_cell", "nk_cell", "monocyte"]


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


# ---------------------------------------------------------------------------
# Part 3 - Statistical analysis: responders vs non-responders
# ---------------------------------------------------------------------------
def responder_frequencies(conn=None, db_path=DB_PATH,
                          condition="melanoma", treatment="miraclib",
                          sample_type="PBMC"):
    """Per-sample population frequencies for the responder comparison subset.

    Returns one row per (sample, population) restricted to the given
    condition / treatment / sample_type and to samples whose subject has a
    known response (yes/no). Columns: sample, subject, response,
    time_from_treatment_start, population, count, total_count, percentage.
    """
    close = conn is None
    conn = conn or connect(db_path)
    try:
        query = """
            SELECT
                s.sample,
                su.subject,
                su.response,
                s.time_from_treatment_start,
                cc.population,
                cc.count,
                SUM(cc.count) OVER (PARTITION BY s.sample)                 AS total_count,
                100.0 * cc.count / SUM(cc.count) OVER (PARTITION BY s.sample) AS percentage
            FROM cell_counts cc
            JOIN samples  s  ON cc.sample  = s.sample
            JOIN subjects su ON s.subject  = su.subject
            WHERE su.condition  = ?
              AND su.treatment  = ?
              AND s.sample_type = ?
              AND su.response IN ('yes', 'no')
            ORDER BY cc.population, su.response, s.sample
        """
        df = pd.read_sql_query(
            query, conn, params=(condition, treatment, sample_type)
        )
    finally:
        if close:
            conn.close()
    return df


def responder_stats(freq_df, alpha=0.05):
    """Compare responders vs non-responders per population.

    For each population, runs a two-sided Mann-Whitney U test (non-parametric:
    relative frequencies are bounded and not guaranteed normal) on the
    responder vs non-responder percentage distributions. p-values are
    corrected for testing 5 populations with the Benjamini-Hochberg FDR
    procedure. Returns a DataFrame sorted by significance.
    """
    rows = []
    for pop in sorted(freq_df["population"].unique()):
        sub = freq_df[freq_df["population"] == pop]
        responders = sub.loc[sub["response"] == "yes", "percentage"]
        non_resp = sub.loc[sub["response"] == "no", "percentage"]
        u_stat, p_val = mannwhitneyu(
            responders, non_resp, alternative="two-sided"
        )
        rows.append({
            "population": pop,
            "n_responder": len(responders),
            "n_non_responder": len(non_resp),
            "median_responder": responders.median(),
            "median_non_responder": non_resp.median(),
            "median_diff": responders.median() - non_resp.median(),
            "u_statistic": u_stat,
            "p_value": p_val,
        })

    stats = pd.DataFrame(rows)
    # Benjamini-Hochberg FDR correction across the 5 populations tested.
    stats["p_value_adj"] = false_discovery_control(stats["p_value"], method="bh")
    stats["significant"] = stats["p_value_adj"] < alpha
    return stats.sort_values("p_value_adj").reset_index(drop=True)


def responder_boxplot(freq_df, outpath="responder_boxplot.png", ax=None):
    """Boxplot of population frequencies, responders vs non-responders.

    Returns the matplotlib Axes. If ``ax`` is None a new figure is created
    and (when ``outpath`` is given) saved to disk.
    """
    import matplotlib
    if ax is None:
        matplotlib.use("Agg")  # headless-safe when called from a script
    import matplotlib.pyplot as plt
    import seaborn as sns

    created = ax is None
    if created:
        fig, ax = plt.subplots(figsize=(10, 6))

    order = sorted(freq_df["population"].unique())
    sns.boxplot(
        data=freq_df, x="population", y="percentage", hue="response",
        order=order, hue_order=["yes", "no"], ax=ax,
    )
    ax.set_xlabel("Immune cell population")
    ax.set_ylabel("Relative frequency (%)")
    ax.set_title("Population relative frequency: responders vs non-responders\n"
                 "(melanoma, miraclib, PBMC)")
    ax.legend(title="Responder")

    if created and outpath:
        ax.figure.tight_layout()
        ax.figure.savefig(outpath, dpi=150)
    return ax


# ---------------------------------------------------------------------------
# Part 4 - Data subset analysis: baseline melanoma/miraclib PBMC samples
# ---------------------------------------------------------------------------
def baseline_subset(conn=None, db_path=DB_PATH, condition="melanoma",
                    treatment="miraclib", sample_type="PBMC", timepoint=0):
    """Return the baseline subset as one row per sample with subject metadata.

    Default filter: melanoma PBMC samples at baseline
    (time_from_treatment_start = 0) from miraclib-treated patients.
    Columns: sample, subject, project, response, sex, age.
    """
    close = conn is None
    conn = conn or connect(db_path)
    try:
        query = """
            SELECT s.sample, su.subject, su.project,
                   su.response, su.sex, su.age
            FROM samples  s
            JOIN subjects su ON s.subject = su.subject
            WHERE su.condition  = ?
              AND su.treatment  = ?
              AND s.sample_type = ?
              AND s.time_from_treatment_start = ?
            ORDER BY s.sample
        """
        df = pd.read_sql_query(
            query, conn, params=(condition, treatment, sample_type, timepoint)
        )
    finally:
        if close:
            conn.close()
    return df


def subset_breakdowns(subset_df):
    """Summarize the baseline subset for Bob's three questions.

    Returns a dict of DataFrames:
        samples_per_project    - sample counts per project
        subjects_by_response   - distinct subject counts (responder/non-resp)
        subjects_by_sex        - distinct subject counts (male/female)
    """
    per_project = (
        subset_df.groupby("project")["sample"].nunique()
        .rename("n_samples").reset_index()
    )
    by_response = (
        subset_df.drop_duplicates("subject").groupby("response")["subject"]
        .nunique().rename("n_subjects").reset_index()
    )
    by_sex = (
        subset_df.drop_duplicates("subject").groupby("sex")["subject"]
        .nunique().rename("n_subjects").reset_index()
    )
    return {
        "samples_per_project": per_project,
        "subjects_by_response": by_response,
        "subjects_by_sex": by_sex,
    }


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

    # -- Part 3: responders vs non-responders (melanoma, miraclib, PBMC) -----
    print("\n" + "=" * 70)
    print("Part 3: responders vs non-responders (melanoma, miraclib, PBMC)")
    print("=" * 70)
    rf = responder_frequencies()
    print(f"Samples in subset: {rf['sample'].nunique():,} "
          f"(responders={rf.loc[rf.response=='yes','sample'].nunique()}, "
          f"non-responders={rf.loc[rf.response=='no','sample'].nunique()})\n")

    stats = responder_stats(rf)
    show = stats.copy()
    for col in ["median_responder", "median_non_responder", "median_diff"]:
        show[col] = show[col].round(2)
    show["p_value"] = show["p_value"].map(lambda v: f"{v:.2e}")
    show["p_value_adj"] = show["p_value_adj"].map(lambda v: f"{v:.2e}")
    print(show.to_string(index=False))

    sig = stats.loc[stats["significant"], "population"].tolist()
    print(f"\nSignificant populations (BH-adjusted p < 0.05): "
          f"{', '.join(sig) if sig else 'none'}")

    fig_path = "responder_boxplot.png"
    responder_boxplot(rf, outpath=fig_path)
    print(f"Boxplot written to {fig_path}")

    # -- Part 4: baseline subset (melanoma, miraclib, PBMC, day 0) ----------
    print("\n" + "=" * 70)
    print("Part 4: baseline subset (melanoma, miraclib, PBMC, time=0)")
    print("=" * 70)
    subset = baseline_subset()
    print(f"Samples: {subset['sample'].nunique():,}  |  "
          f"Subjects: {subset['subject'].nunique():,}\n")
    breakdowns = subset_breakdowns(subset)
    print("Samples per project:")
    print(breakdowns["samples_per_project"].to_string(index=False))
    print("\nSubjects by response:")
    print(breakdowns["subjects_by_response"].to_string(index=False))
    print("\nSubjects by sex:")
    print(breakdowns["subjects_by_sex"].to_string(index=False))
