#!/usr/bin/env python3
"""Interactive dashboard for the cell-count analysis (Loblaw Bio).

Run with:

    python app.py

then open the printed URL (default http://127.0.0.1:8050) in a browser.

Three tabs map to the analysis tasks:
  - Overview (Part 2): per-sample relative frequency summary table.
  - Responders (Part 3): responder vs non-responder boxplot + statistics,
    with configurable filters and an optional single-timepoint view.
  - Baseline subset (Part 4): filterable subset breakdowns by project,
    response, and sex.

All data comes from analysis.py, which queries cell-count.db.
"""

import os

import pandas as pd
import plotly.express as px
import plotly.graph_objects as go
from dash import Dash, dash_table, dcc, html, Input, Output

import analysis

DB_PATH = "cell-count.db"
POP_ORDER = ["b_cell", "cd4_t_cell", "cd8_t_cell", "monocyte", "nk_cell"]


def _require_db():
    if not os.path.exists(DB_PATH):
        raise SystemExit(
            f"{DB_PATH} not found. Run 'python load_data.py' first."
        )


_require_db()
OPTIONS = analysis.filter_options()


# ---------------------------------------------------------------------------
# Reusable UI helpers
# ---------------------------------------------------------------------------
def dropdown(id_, values, value, clearable=False):
    return dcc.Dropdown(
        id=id_,
        options=[{"label": str(v), "value": v} for v in values],
        value=value,
        clearable=clearable,
        style={"width": "180px"},
    )


def data_table(id_, page_size=15):
    return dash_table.DataTable(
        id=id_,
        page_size=page_size,
        sort_action="native",
        filter_action="native",
        style_table={"overflowX": "auto"},
        style_cell={"fontFamily": "sans-serif", "padding": "6px",
                    "textAlign": "left"},
        style_header={"fontWeight": "bold", "backgroundColor": "#f0f2f6"},
    )


def _df_to_records(df, round_cols=None, fmt=None):
    out = df.copy()
    if round_cols:
        for c in round_cols:
            out[c] = out[c].round(2)
    if fmt:
        for c, f in fmt.items():
            out[c] = out[c].map(f)
    return out.to_dict("records"), [{"name": c, "id": c} for c in out.columns]


# ---------------------------------------------------------------------------
# App layout
# ---------------------------------------------------------------------------
app = Dash(__name__, title="Loblaw Bio — Cell Count Dashboard")
server = app.server  # for optional WSGI deployment

filter_row_style = {"display": "flex", "gap": "24px", "flexWrap": "wrap",
                    "alignItems": "flex-end", "marginBottom": "16px"}

# Each Part 4 chart column shares the row equally and is allowed to shrink
# (minWidth 0) so all three stay side by side instead of wrapping/overflowing.
_chart_col_style = {"flex": "1 1 0", "minWidth": 0, "overflow": "hidden"}


def labeled(label, component):
    return html.Div([html.Label(label, style={"fontWeight": "bold",
                                               "display": "block"}),
                     component])


# ----- Tab 1: Overview (Part 2) --------------------------------------------
overview_tab = html.Div([
    html.H3("Per-sample cell population relative frequencies"),
    html.P("Each row is one population in one sample: total_count is the sum "
           "across all five populations, and percentage is the population's "
           "share of that total."),
    html.Div(labeled("Filter by sample id (optional)",
                     dcc.Input(id="ov-sample", type="text", value="",
                               placeholder="e.g. sample00000",
                               style={"width": "220px"}))),
    html.Br(),
    html.Div(id="ov-count", style={"marginBottom": "8px",
                                   "fontStyle": "italic"}),
    data_table("ov-table"),
], style={"padding": "16px"})


# ----- Tab 2: Responders (Part 3) ------------------------------------------
responders_tab = html.Div([
    html.H3("Responders vs non-responders"),
    html.P("Compares population relative frequencies between responders and "
           "non-responders. Significance: two-sided Mann-Whitney U per "
           "population, Benjamini-Hochberg FDR corrected (alpha = 0.05)."),
    html.Div([
        labeled("Condition", dropdown("rp-condition", OPTIONS["condition"],
                                      "melanoma")),
        labeled("Treatment", dropdown("rp-treatment", OPTIONS["treatment"],
                                      "miraclib")),
        labeled("Sample type", dropdown("rp-sample-type",
                                        OPTIONS["sample_type"], "PBMC")),
        labeled("Timepoint", dropdown(
            "rp-timepoint", ["All"] + OPTIONS["timepoint"], "All")),
    ], style=filter_row_style),
    html.Div(id="rp-summary", style={"marginBottom": "8px",
                                     "fontStyle": "italic"}),
    dcc.Graph(id="rp-boxplot"),
    html.H4("Per-population statistics"),
    data_table("rp-stats", page_size=10),
], style={"padding": "16px"})


# ----- Tab 3: Baseline subset (Part 4) -------------------------------------
subset_tab = html.Div([
    html.H3("Data subset breakdowns"),
    html.P("Defaults to the Part 4 question: baseline (day 0) melanoma PBMC "
           "samples from miraclib-treated patients. Adjust the filters to "
           "explore other subsets."),
    html.Div([
        labeled("Condition", dropdown("sb-condition", OPTIONS["condition"],
                                      "melanoma")),
        labeled("Treatment", dropdown("sb-treatment", OPTIONS["treatment"],
                                      "miraclib")),
        labeled("Sample type", dropdown("sb-sample-type",
                                        OPTIONS["sample_type"], "PBMC")),
        labeled("Timepoint", dropdown("sb-timepoint", OPTIONS["timepoint"],
                                      0)),
    ], style=filter_row_style),
    html.Div(id="sb-summary", style={"fontWeight": "bold",
                                     "marginBottom": "12px"}),
    html.Div([
        html.Div([html.H4("Samples per project"),
                  dcc.Graph(id="sb-project", style={"height": "320px"},
                            config={"responsive": True})],
                 style=_chart_col_style),
        html.Div([html.H4("Subjects by response"),
                  dcc.Graph(id="sb-response", style={"height": "320px"},
                            config={"responsive": True})],
                 style=_chart_col_style),
        html.Div([html.H4("Subjects by sex"),
                  dcc.Graph(id="sb-sex", style={"height": "320px"},
                            config={"responsive": True})],
                 style=_chart_col_style),
    ], style={"display": "flex", "gap": "16px", "flexWrap": "nowrap",
              "width": "100%"}),
], style={"padding": "16px"})


app.layout = html.Div([
    html.H1("Loblaw Bio — Immune Cell Count Dashboard"),
    html.P("Exploring how miraclib affects immune cell populations across "
           "patient samples."),
    dcc.Tabs([
        dcc.Tab(label="Overview (Part 2)", children=overview_tab),
        dcc.Tab(label="Responders (Part 3)", children=responders_tab),
        dcc.Tab(label="Baseline subset (Part 4)", children=subset_tab),
    ]),
], style={"maxWidth": "1200px", "margin": "0 auto",
          "fontFamily": "sans-serif"})


# ---------------------------------------------------------------------------
# Callbacks
# ---------------------------------------------------------------------------
@app.callback(
    Output("ov-table", "data"),
    Output("ov-table", "columns"),
    Output("ov-count", "children"),
    Input("ov-sample", "value"),
)
def update_overview(sample_query):
    freq = analysis.cell_frequencies(db_path=DB_PATH)
    if sample_query:
        freq = freq[freq["sample"].str.contains(sample_query.strip(),
                                                case=False, na=False)]
    data, columns = _df_to_records(freq, round_cols=["percentage"])
    msg = (f"{len(freq):,} rows across {freq['sample'].nunique():,} samples"
           if len(freq) else "No samples match that filter.")
    return data, columns, msg


@app.callback(
    Output("rp-boxplot", "figure"),
    Output("rp-stats", "data"),
    Output("rp-stats", "columns"),
    Output("rp-summary", "children"),
    Input("rp-condition", "value"),
    Input("rp-treatment", "value"),
    Input("rp-sample-type", "value"),
    Input("rp-timepoint", "value"),
)
def update_responders(condition, treatment, sample_type, timepoint):
    tp = None if timepoint == "All" else timepoint
    rf = analysis.responder_frequencies(
        db_path=DB_PATH, condition=condition, treatment=treatment,
        sample_type=sample_type, timepoint=tp,
    )

    n_yes = rf.loc[rf["response"] == "yes", "sample"].nunique()
    n_no = rf.loc[rf["response"] == "no", "sample"].nunique()

    # Need both groups (and >1 value) to draw a meaningful comparison.
    if n_yes == 0 or n_no == 0:
        empty = go.Figure().update_layout(
            title="No responder/non-responder samples for this filter.")
        return empty, [], [], (
            f"Responders: {n_yes} samples | Non-responders: {n_no} samples "
            "— need both groups to compare.")

    tp_label = "all timepoints" if tp is None else f"day {tp}"
    fig = px.box(
        rf, x="population", y="percentage", color="response",
        category_orders={"population": POP_ORDER,
                         "response": ["yes", "no"]},
        color_discrete_map={"yes": "#1f77b4", "no": "#ff7f0e"},
        labels={"percentage": "Relative frequency (%)",
                "population": "Immune cell population",
                "response": "Responder"},
        title=f"{condition} / {treatment} / {sample_type} ({tp_label})",
    )
    fig.update_layout(boxmode="group")

    stats = analysis.responder_stats(rf)
    data, columns = _df_to_records(
        stats,
        round_cols=["median_responder", "median_non_responder",
                    "median_diff"],
        fmt={"p_value": lambda v: f"{v:.2e}",
             "p_value_adj": lambda v: f"{v:.2e}",
             "u_statistic": lambda v: f"{v:.0f}"},
    )
    sig = stats.loc[stats["significant"], "population"].tolist()
    summary = (
        f"Responders: {n_yes} samples | Non-responders: {n_no} samples. "
        f"Significant populations (adj p < 0.05): "
        f"{', '.join(sig) if sig else 'none'}.")
    return fig, data, columns, summary


def _bar(df, x, y, color):
    fig = px.bar(df, x=x, y=y, text=y)
    fig.update_traces(marker_color=color, textposition="outside")
    # autosize lets the figure fill the fixed-height, shrinkable column rather
    # than imposing Plotly's default ~450px height.
    fig.update_layout(autosize=True, height=None,
                      margin={"t": 20, "b": 30, "l": 40, "r": 10},
                      yaxis_title=y, xaxis_title="", showlegend=False)
    return fig


@app.callback(
    Output("sb-summary", "children"),
    Output("sb-project", "figure"),
    Output("sb-response", "figure"),
    Output("sb-sex", "figure"),
    Input("sb-condition", "value"),
    Input("sb-treatment", "value"),
    Input("sb-sample-type", "value"),
    Input("sb-timepoint", "value"),
)
def update_subset(condition, treatment, sample_type, timepoint):
    subset = analysis.baseline_subset(
        db_path=DB_PATH, condition=condition, treatment=treatment,
        sample_type=sample_type, timepoint=timepoint,
    )
    summary = (f"{subset['sample'].nunique():,} samples from "
               f"{subset['subject'].nunique():,} subjects "
               f"({condition} / {treatment} / {sample_type} / day {timepoint})")

    if subset.empty:
        blank = go.Figure().update_layout(title="No samples")
        return summary + " — no samples match this filter.", blank, blank, blank

    b = analysis.subset_breakdowns(subset)
    return (
        summary,
        _bar(b["samples_per_project"], "project", "n_samples", "#4c78a8"),
        _bar(b["subjects_by_response"], "response", "n_subjects", "#54a24b"),
        _bar(b["subjects_by_sex"], "sex", "n_subjects", "#e45756"),
    )


if __name__ == "__main__":
    app.run(debug=True)
