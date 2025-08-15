import pandas as pd
from scipy.stats import chi2_contingency
import numpy as np
from statsmodels.stats.multitest import multipletests
import plotly.graph_objects as go
import plotly.express as px

# Function to retrieve metrics (total counts, % of transcript with at least 1 count)
def get_metrics(input_df: pd.DataFrame, txid_results: pd.DataFrame = None) -> pd.DataFrame:
    """
    Function that calculates, for each Run, the total counts and the percentage of the transcript
    having at least one count.

    :param all_studies_df: Dataframe with all the raw counts
    :param txid_results: Optional. Dataframe with the results of the main analysis. If none is given, metrics
    will be calculated for all studies, not just results.
    :return: Dataframe with the total counts and % for each individual run (SRR)
    """
    df = input_df.copy()

    if txid_results is not None:
        good_runs = txid_results.index.to_list()
        df = input_df[good_runs].copy()

    # Code to calculate the % covered
    df2 = df.copy()
    df2 = (df2 != 0).astype(int)
    total_percentage = df2.sum(axis=0)/df2.shape[0]
    total_percentage.name = "%_covered"

    # Code to Retrieve Total Counts
    totals = df.sum(axis=0)
    totals.name = "tot_counts"

    # Set up results dataframe
    metrics_df = pd.concat([totals, total_percentage], axis=1)
    metrics_df["Run"] = metrics_df.index
    new_order = ["Run", "tot_counts","%_covered"]  # your desired order
    metrics_df = metrics_df[new_order]

    return metrics_df


# Function to do the breakdown from the full transcript
def get_breakdowns_full_tx(all_studies_df: pd.DataFrame,
                           utr_5: int,
                           cds: int,
                           txid_results: pd.DataFrame = None):
    """
    :param all_studies_df: Dataframe with all the raw counts
    :param utr_5: length of 5' utr
    :param cds: length of CDS
    :param txid_results: Optional. Dataframe with the results of the main analysis. If none is given, metrics
    will be calculated for all studies, not just results.
    :return: Dataframe with total counts for each region (utrs, cds)
    """

    df = all_studies_df.copy()
    if txid_results is not None:
        good_runs = txid_results.index.to_list()
        df = all_studies_df[good_runs].copy()

    # Compute global aggregate
    df["Global Aggregate"] = df.sum(axis=1)
    # Compute the three totals by region...
    leader_totals = df.iloc[:utr_5].sum(axis=0)
    leader_totals.name = "UTR_5"
    cds_totals = df.iloc[utr_5:utr_5 + cds].sum(axis=0)
    cds_totals.name = "CDS"
    trailer_totals = df.iloc[utr_5 + cds:].sum(axis=0)
    trailer_totals.name = "UTR_3"

    breakdown = pd.concat([leader_totals, cds_totals, trailer_totals], axis=1)
    raw_counts_breakdown = breakdown.copy()

    breakdown["tot_counts"] = breakdown["UTR_5"] + breakdown["CDS"] + breakdown["UTR_3"]
    breakdown["%_UTR_5"] = breakdown["UTR_5"]/breakdown["tot_counts"]
    breakdown["%_CDS"] = breakdown["CDS"]/breakdown["tot_counts"]
    breakdown["%_UTR_3"] = breakdown["UTR_3"]/breakdown["tot_counts"]

    breakdown["Run"] = breakdown.index
    # Add part that computes totals and percentage in each section and we are golden
    cols_to_keep = ["Run", "tot_counts", "%_UTR_5", "%_CDS", "%_UTR_3" ]
    totals_and_percents_breakdown = breakdown[cols_to_keep].copy()

    return raw_counts_breakdown, totals_and_percents_breakdown


# Function to do breakdowns across (user-selected) regions from the raw counts.
# Also does the percentage breakdowns?
def get_breakdowns_slices(all_studies_df: pd.DataFrame,
                          coordinates: dict,
                          txid_results: pd.DataFrame = None):
    """
    :param all_studies_df: Dataframe with all the raw counts
    :param coordinates: Dictionary containing info on how to make breakdowns. a dict of tuples, with the slices
    used in the analysis ({x_slice: (,), y_slice: (,), n_slice: (,)})
    :param txid_results: Optional. Dataframe with the results of the main analysis. If none is given, metrics
    will be calculated for all studies, not just results.
    :return: a dataframe with the total counts per each region (the breakdowns) across all SRRs.
    """

    df = all_studies_df.copy()

    if txid_results is not None:
        good_runs = txid_results.index.to_list()
        df = all_studies_df[good_runs].copy()

    # Slice the dataframe, using the coordinates, in three subsets.
    slice_x = df.iloc[coordinates["x_slice"][0]:coordinates["x_slice"][1]].copy()
    slice_y = df.iloc[coordinates["y_slice"][0]:coordinates["y_slice"][1]].copy()
    slice_n = df.iloc[coordinates["n_slice"][0]:coordinates["n_slice"][1]].copy()

    # Compute metrics for each (total of each region, Global Aggregate, percentage of total counts)
    x_totals = slice_x.sum(axis=0).rename("First Region")
    x_totals["Global Aggregate"] = x_totals.sum() # Does this aggregate like I want?

    y_totals = slice_y.sum(axis=0).rename("Second Region")
    y_totals["Global Aggregate"] = y_totals.sum()

    n_totals = slice_n.sum(axis=0).rename("Noise Region")
    n_totals["Global Aggregate"] = n_totals.sum()

    # Merge again into a single one
    breakdown = pd.concat([x_totals, y_totals, n_totals], axis = 1)
    breakdown["tot_counts"] = breakdown["First Region"] + breakdown["Second Region"] + breakdown["Noise Region"]

    # We just add the value, we will divide later to obtain the percentage
    breakdown["% First Region"] = breakdown["First Region"]/breakdown["tot_counts"]
    breakdown["% Second Region"] = breakdown["Second Region"]/breakdown["tot_counts"]
    breakdown["% Noise Region"] = breakdown["Noise Region"]/breakdown["tot_counts"]

    # Divide the breakdowns in two dataframes, one will be used in X-squared, the other in results plotting
    cols_for_X_test = ["First Region", "Second Region", "Noise Region"]
    raw_counts_breakdowns = breakdown[cols_for_X_test].copy()

    cols_for_plots = ["tot_counts", "% First Region","% Second Region", "% Noise Region"]
    totals_and_percents_breakdowns = breakdown[cols_for_plots].copy()
    totals_and_percents_breakdowns["Run"] = totals_and_percents_breakdowns.index

    # Extract global aggregate
    global_aggregate = breakdown.loc[["Global Aggregate"]].copy()
    global_aggregate =  global_aggregate.drop(columns=["First Region", "Second Region", "Noise Region"])
    global_aggregate["Run"] = "Global Aggregate"

    return raw_counts_breakdowns, totals_and_percents_breakdowns, global_aggregate


# Function to perform the pair-wise X-squared test, with also Cramer's V, Bonferroni and BH Corrections.
def perform_pairwise_chi(df2, baseline_idx):
    """
    Perform pairwise Chi-squared tests of all rows in df2 versus the baseline row.
    Returns a DataFrame with chi2, p-value, dof, and corrected p-values.
    """

    results = []

    for idx, row in df2.iterrows():
        if idx == baseline_idx:
            continue

        table = pd.DataFrame([df2.loc[baseline_idx], row])
        chi2, p, dof, expected = chi2_contingency(table)
        n = table.to_numpy().sum()
        r, k = table.shape
        cramers_v = np.sqrt(chi2 / (n * (min(k - 1, r - 1))))
        results.append({
            "Run": idx,
            "chi2": chi2,
            "p-value": p,
            "dof": dof,
            "cramers_v": cramers_v
        })

    results_df = pd.DataFrame(results)

    if not results_df.empty:
        # Bonferroni correction
        bonf_corrected = multipletests(results_df['p-value'], method='bonferroni')[1]
        # Benjamini-Hochberg FDR correction
        fdr_corrected = multipletests(results_df['p-value'], method='fdr_bh')[1]

        results_df['bonferroni_p'] = bonf_corrected
        results_df['fdr_p'] = fdr_corrected

    return results_df

# Function to prepare results for plotter
def prepare_chisquared_results(results_df, percent_breakdowns, global_aggregate_values, correction = "bonferroni_p"):
    # Merge dataframes
    df = pd.merge(results_df, percent_breakdowns, on="Run")
    # Remove n.s. results
    cleaned_df = df[df[correction] < 0.05].copy()

    # Reintroduce the Global Aggregate
    global_aggregate_values["cramers_v"] = 0
    global_aggregate_values[correction] = 1

    # Concatenate to reintroduce the Global Aggregate after the analysis is done
    df = pd.concat([cleaned_df, global_aggregate_values], ignore_index=True)
    df = df.where(pd.notnull(df), None)

    # Order by tot_counts
    df = df.sort_values(by="tot_counts", ascending=True)

    return df


# Function to make plots

def better_plotter(df: pd.DataFrame, region_cols):
    x_labels = region_cols
    y_labels = df["Run"].astype(str)
    z = df[region_cols].values  # n_rows x 3

    # Hover text: include both region value and p-value for the row
    hovertext = np.array([
        [f"{row}<br>{col}: {val:.2f}% of Translational Activity <br>Adjusted (Bonf) p-value: {pval:.3g}, Total counts: {tot_counts}"
         for col, val in zip(region_cols, row_vals)]
        for row, row_vals, pval, tot_counts in zip(y_labels, z, df['bonferroni_p'], df['tot_counts'])
    ])

    fig = go.Figure()

    # Heatmap
    fig.add_trace(go.Heatmap(
        z=z,
        x=x_labels,
        y=y_labels,
        colorscale='Viridis',
        zmin=0, zmax=1,
        hoverinfo='text',
        text=hovertext,
        showscale=True,
        colorbar=dict(
            title='% Region',
            thickness = 12,
            len = 0.1,
            y = 0.9)
    ))
    # === Cramér’s V text "column" on the right ===
    if 'cramers_v' in df.columns:
        v_text = [f"V={v:.2f}" if pd.notna(v) else "V=NA" for v in df['cramers_v']]
        v_hover = [f"{run}<br>Cramér’s V: {v:.3f}" if pd.notna(v) else f"{run}<br>Cramér’s V: NA"
                   for run, v in zip(y_labels, df['cramers_v'])]

        fig.add_trace(go.Scatter(
            x=['cramers_v'] * len(df),          # new categorical column to the right
            y=y_labels,
            mode='text',
            text=v_text,
            hoverinfo='text',
            hovertext=v_hover,
            textfont=dict(size=12),
            showlegend=False
        ))

        # Ensure the x-axis shows columns in desired order
        fig.update_xaxes(categoryorder='array',
                         categoryarray=x_labels + ['cramers_v'])




    fig.update_layout(
        title="Heatmap of Regional Translational Acitivity (%)<br> <br>",
        xaxis_title="",
        yaxis_title="Run",
        width=420,
        height=30 * len(df),
        margin=dict(l=60, r=40, t=60, b=40)
    )

    fig.update_xaxes(side='top')

    return fig

